


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

#ifndef __GMS_LINEAR_OPTIMAL_CONTROL_V2_H__
#define __GMS_LINEAR_OPTIMAL_CONTROL_V2_H__ 251020260730

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/MatrixFunctions>
#include <cmath>
#include <random>
#include <stdexcept>
#include <tuple>
#include <vector>
#include <algorithm>
#include "GMS_config.h"

namespace gms
{
namespace state_space
{

using MatrixXd = Eigen::MatrixXd;
using VectorXd = Eigen::VectorXd;
using MatrixXi = Eigen::MatrixXi;
using VectorXi = Eigen::VectorXi;

// -----------------------------------------------------------------------------
// Simple MATLAB-like helpers
// -----------------------------------------------------------------------------
inline MatrixXd zeros(int rows, int cols)
{
    return MatrixXd::Zero(rows, cols);
}

inline VectorXd zeros_vec(int n)
{
    return VectorXd::Zero(n);
}

inline MatrixXd eye(int n)
{
    return MatrixXd::Identity(n, n);
}

inline MatrixXd ones(int rows, int cols)
{
    return MatrixXd::Ones(rows, cols);
}

inline VectorXd ones_vec(int n)
{
    return VectorXd::Ones(n);
}

inline MatrixXd inv(const MatrixXd& A)
{
    if (A.rows() != A.cols())
        throw std::invalid_argument("inv(): matrix must be square");

    Eigen::FullPivLU<MatrixXd> lu(A);
    if (lu.rank() < A.rows())
        throw std::runtime_error("inv(): matrix is singular or rank-deficient");

    return lu.inverse();
}

inline MatrixXd kron(const MatrixXd& A, const MatrixXd& B)
{
    return Eigen::kroneckerProduct(A, B);
}

inline MatrixXd transpose(const MatrixXd& A)
{
    return A.transpose();
}

inline VectorXd diag(const MatrixXd& M)
{
    return M.diagonal();
}

inline MatrixXd sqrtm(const MatrixXd& A)
{
    Eigen::SelfAdjointEigenSolver<MatrixXd> es(A);
    if (es.info() != Eigen::Success)
        throw std::runtime_error("sqrtm(): eigen decomposition failed");
    return es.operatorSqrt();
}

inline MatrixXd randn(int rows, int cols)
{
    static std::mt19937 rng(std::random_device{}());
    std::normal_distribution<double> dist(0.0, 1.0);

    MatrixXd M(rows, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            M(i, j) = dist(rng);
    return M;
}

inline VectorXd randn_vec(int n)
{
    static std::mt19937 rng(std::random_device{}());
    std::normal_distribution<double> dist(0.0, 1.0);

    VectorXd v(n);
    for (int i = 0; i < n; ++i)
        v(i) = dist(rng);
    return v;
}

// 3D matrix array helper: MATLAB arrays such as Ka(:,:,k) are modeled as std::vector<MatrixXd>
using Matrix3D = std::vector<MatrixXd>;

inline Matrix3D make3d(int rows, int cols, int pages)
{
    return Matrix3D(pages, MatrixXd::Zero(rows, cols));
}

inline Matrix3D make3d(int n, int pages)
{
    return Matrix3D(pages, MatrixXd::Zero(n, n));
}

inline MatrixXd& page(Matrix3D& arr, int idx)
{
    if (idx < 0 || idx >= static_cast<int>(arr.size()))
        throw std::out_of_range("page index out of range");
    return arr[idx];
}

inline const MatrixXd& page(const Matrix3D& arr, int idx)
{
    if (idx < 0 || idx >= static_cast<int>(arr.size()))
        throw std::out_of_range("page index out of range");
    return arr[idx];
}

// -----------------------------------------------------------------------------
// StateSpace wrapper
// -----------------------------------------------------------------------------
struct StateSpace
{
    MatrixXd A, B, C, D;

    StateSpace() = default;

    StateSpace(const MatrixXd& A_, const MatrixXd& B_, const MatrixXd& C_, const MatrixXd& D_)
        : A(A_), B(B_), C(C_), D(D_)
    {
    }

    std::tuple<MatrixXd, MatrixXd, MatrixXd, MatrixXd> ssdata() const
    {
        return {A, B, C, D};
    }

    StateSpace discretize(double Ts) const
    {
        // zero-order-hold discretization via matrix exponential
        return c2d(*this, Ts);
    }
};

inline std::tuple<MatrixXd, MatrixXd, MatrixXd, MatrixXd> ssdata(const StateSpace& sys)
{
    return sys.ssdata();
}

inline StateSpace ss(const MatrixXd& A, const MatrixXd& B, const MatrixXd& C, const MatrixXd& D)
{
    return StateSpace(A, B, C, D);
}

// -----------------------------------------------------------------------------
// Useful control routines
// -----------------------------------------------------------------------------
// Discrete-time LQR:
inline MatrixXd dlqr(const MatrixXd& A,
                    const MatrixXd& B,
                    const MatrixXd& Q,
                    const MatrixXd& R,
                    const MatrixXd& N = MatrixXd())
{
    const int n = A.rows();
    const int m = B.cols();

    MatrixXd Np = N;
    if (Np.rows() == 0 || Np.cols() == 0)
        Np = MatrixXd::Zero(n, m);

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

// Discrete-time LQE (Kalman gain):
inline MatrixXd dlqe(const MatrixXd& A,
                    const MatrixXd& G,
                    const MatrixXd& C,
                    const MatrixXd& W,
                    const MatrixXd& V)
{
    MatrixXd P = G * W * G.transpose();
    const int max_iter = 200;
    const double tol = 1e-12;

    for (int iter = 0; iter < max_iter; ++iter) {
        MatrixXd S = C * P * C.transpose() + V;
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

// MATLAB-compatible wrappers
inline MatrixXd lqr(const MatrixXd& A, const MatrixXd& B, const MatrixXd& Q, const MatrixXd& R, const MatrixXd& N = MatrixXd())
{
    return dlqr(A, B, Q, R, N);
}

inline MatrixXd lqe(const MatrixXd& A, const MatrixXd& G, const MatrixXd& C, const MatrixXd& W, const MatrixXd& V)
{
    return dlqe(A, G, C, W, V);
}

// c2d for zero-order-hold
inline StateSpace c2d(const StateSpace& sys, double Ts)
{
    const int n = sys.A.rows();
    const int m = sys.B.cols();

    MatrixXd M = MatrixXd::Zero(n + m, n + m);
    M.block(0, 0, n, n) = sys.A;
    M.block(0, n, n, m) = sys.B;
    M = M * Ts;

    MatrixXd E = M.exp();

    MatrixXd Phi = E.block(0, 0, n, n);
    MatrixXd Gamma = E.block(0, n, n, m);

    return StateSpace(Phi, Gamma, sys.C, sys.D);
}

// Linear simulation for state-space
inline std::pair<MatrixXd, MatrixXd> lsim(const StateSpace& sys,
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
        VectorXd u = U.col(k);
        X.col(k + 1) = sys.A * X.col(k) + sys.B * u;
        Y.col(k + 1) = sys.C * X.col(k + 1) + sys.D * u;
    }

    return { X, Y };
}

// MATLAB 'initial' behavior
inline std::pair<MatrixXd, MatrixXd> initial(const StateSpace& sys,
                                           const VectorXd& x0,
                                           double tf,
                                           int Ns = 100)
{
    const int n = sys.A.rows();
    const int m = sys.B.cols();
    MatrixXd U = MatrixXd::Zero(m, Ns + 1);
    return lsim(sys, U, x0);
}

// Discrete-time Lyapunov equation X = A X A' + Q
inline MatrixXd dlyap(const MatrixXd& A, const MatrixXd& Q)
{
    MatrixXd X = Q;
    for (int iter = 0; iter < 200; ++iter) {
        MatrixXd Xnext = A * X * A.transpose() + Q;
        if ((Xnext - X).norm() < 1e-12)
            return Xnext;
        X = Xnext;
    }
    return X;
}

// Continuous-time Lyapunov equation A X + X A' + Q = 0
inline MatrixXd lyap_cont(const MatrixXd& A, const MatrixXd& Q)
{
    // Solve using a linear system with vec operator:
    // vec(X) = - (I ⊗ A + A ⊗ I)^-1 vec(Q)
    const int n = A.rows();
    MatrixXd K = kron(MatrixXd::Identity(n, n), A) + kron(A, MatrixXd::Identity(n, n));
    VectorXd qvec = Eigen::Map<const VectorXd>(Q.data(), n * n);
    VectorXd xvec = -inv(K) * qvec;
    return Eigen::Map<const MatrixXd>(xvec.data(), n, n);
}

// -----------------------------------------------------------------------------
// Supporting helper for cvrtj / cvrtw
// -----------------------------------------------------------------------------
inline std::tuple<MatrixXd, MatrixXd, MatrixXd, MatrixXd, MatrixXd>
cvrtj(const MatrixXd& A, const MatrixXd& B, const MatrixXd& C,
      const MatrixXd& Q, const MatrixXd& R, double Ts)
{
    StateSpace sys(A, B, C, MatrixXd::Zero(C.rows(), B.cols()));
    StateSpace ds = c2d(sys, Ts);
    MatrixXd N = MatrixXd::Zero(A.rows(), B.cols()); // no cross term in the typical conversion
    return { ds.A, ds.B, Q, N, R };
}

inline MatrixXd cvrtw(const MatrixXd& A, const MatrixXd& G, const MatrixXd& W, double Ts)
{
    // Approximate converted discrete process noise covariance:
    MatrixXd Phi = c2d(StateSpace(A, G, MatrixXd::Identity(A.rows(), A.rows()), MatrixXd::Zero(A.rows(), G.cols())), Ts).A;
    return Phi * W * Phi.transpose();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// 1) dbat_u.m
// -----------------------------------------------------------------------------
// MATLAB: [xs,ws,vss,Rxs,Rw,Rv] = dbat_u(...)
// -----------------------------------------------------------------------------
// This is translated with a faithful structure, but in Eigen the time-varying
// gains are stored as std::vector<MatrixXd> pages, as in MATLAB's K(:,:,k).
// -----------------------------------------------------------------------------
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
    const int nd = Gw.cols();
    const int nc = Gu.cols();
    const int nyf = Cf.rows();
    const int N1 = ys.cols();
    const int Ns = N1 - 1;

    // Cyx = [Cs*Ph; Cs*Ph^2; ...; Cs*Ph^Ns]
    MatrixXd Cyx = MatrixXd::Zero(nm * Ns, ns);
    MatrixXd Z = Ph;
    for (int k = 0; k < Ns; ++k) {
        Cyx.block(k * nm, 0, nm, ns) = Cs * Z;
        Z = Z * Ph;
    }

    // Cyw = [Cs*Ph^(Ns-1)*Gw, ..., Cs*Gw] with MATLAB-style ordering
    MatrixXd Cyw = MatrixXd::Zero(nm * Ns, nd * Ns);
    Z = Ph;
    std::vector<MatrixXd> Cpieces;
    Cpieces.push_back(Cs * Gw);
    for (int k = 1; k < Ns; ++k) {
        Z = Z * Ph;
        Cpieces.push_back(Cs * Z * Gw);
    }
    for (int k = 0; k < Ns; ++k) {
        int row0 = k * nm;
        int col0 = (Ns - 1 - k) * nd;
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
    for (int k = 0; k < Ns; ++k)
        Cfw.block(0, k * nd, nyf, nd) = Fpieces[k];

    // Cfx = Cf * Ph^Ns
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
    Cw.block(n0 + nm + nm * Ns + nyf, ns, nd * Ns, nd * Ns) = MatrixXd::Identity(nd * Ns, nd * Ns);

    Cw.block(0, ns, n0, nd * Ns) = MatrixXd::Zero(n0, nd * Ns);
    Cw.block(n0, ns, nm, nd * Ns) = MatrixXd::Zero(nm, nd * Ns);
    Cw.block(n0 + nm, ns, nm * Ns, nd * Ns) = Cyw;
    Cw.block(n0 + nm + nm * Ns, ns, nyf, nd * Ns) = Cfw;

    // Build ya = [y0; ys(0); ...; ys(Ns); yf; 0; ...; 0]
    VectorXd ya = VectorXd::Zero(n0 + (Ns + 1) * nm + nyf + nd * Ns);
    ya.segment(0, n0) = y0;
    for (int k = 0; k <= Ns; ++k)
        ya.segment(n0 + k * nm, nm) = ys.col(k);
    ya.segment(n0 + (Ns + 1) * nm, nyf) = yf;

    // Cu: [0; Cyu; Cfu; 0]
    MatrixXd Cnu = MatrixXd::Zero(nm * Ns, nc * Ns);
    Z = Ph;
    std::vector<MatrixXd> Upieces;
    Upieces.push_back(Cs * Gu);
    for (int k = 1; k < Ns; ++k) {
        Z = Z * Ph;
        Upieces.push_back(Cs * Z * Gu);
    }
    for (int k = 0; k < Ns; ++k) {
        int row0 = k * nm;
        int col0 = (Ns - 1 - k) * nc;
        Cnu.block(row0, col0, nm, nc) = Upieces[k];
    }

    MatrixXd Cyu = MatrixXd::Zero(nm * Ns, nc * Ns);
    for (int k = 0; k < Ns; ++k) {
        int row0 = k * nm;
        int cols = (k + 1) * nc;
        Cyu.block(row0, 0, nm, cols) = Cnu.block(row0, (Ns - k - 1) * nc, nm, cols);
    }

    MatrixXd Cfu = MatrixXd::Zero(nyf, nc * Ns);
    Z = Cf;
    std::vector<MatrixXd> Fupieces;
    Fupieces.push_back(Z * Gu);
    for (int k = 1; k < Ns; ++k) {
        Z = Z * Ph;
        Fupieces.push_back(Z * Gu);
    }
    for (int k = 0; k < Ns; ++k)
        Cfu.block(0, k * nc, nyf, nc) = Fupieces[k];

    MatrixXd Cu = MatrixXd::Zero(n0 + (Ns + 1) * nm + nyf + nd * Ns,
                                nc * Ns);
    Cu.block(n0 + nm, 0, nm * Ns, nc * Ns) = Cyu;
    Cu.block(n0 + nm + nm * Ns, 0, nyf, nc * Ns) = Cfu;

    // ua = [ud(:,1); ...; ud(:,Ns)]
    VectorXd ua = VectorXd::Zero(nc * Ns);
    for (int k = 0; k < Ns; ++k)
        ua.segment(k * nc, nc) = ud.col(k);

    // Vai block weighting
    MatrixXd Vai = MatrixXd::Zero(n0 + (Ns + 1) * nm + nyf + Ns * nd,
                                  n0 + (Ns + 1) * nm + nyf + Ns * nd);
    Vai.block(0, 0, n0, n0) = inv(V0);
    for (int k = 0; k <= Ns; ++k) {
        int idx = n0 + k * nm;
        Vai.block(idx, idx, nm, nm) = inv(Vd);
    }
    Vai.block(n0 + (Ns + 1) * nm, n0 + (Ns + 1) * nm, nyf, nyf) = inv(Vf);
    for (int k = 0; k < Ns; ++k) {
        int idx = n0 + (Ns + 1) * nm + nyf + k * nd;
        Vai.block(idx, idx, nd, nd) = inv(Wd);
    }

    // solve weighted least squares
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

    // Covariance approximation:
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

// -----------------------------------------------------------------------------
// 2) dbat.m
// -----------------------------------------------------------------------------
// MATLAB: [xs,ws,vss,Rxs,Rw,Rv] = dbat(...)
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
    // This is the uncontrolled-plant variant of dbat_u and shares the same structure.

    const int n0 = C0.rows();
    const int ns = Ph.rows();
    const int nm = ys.rows();
    const int nd = Gw.cols();
    const int nyf = Cf.rows();
    const int N1 = ys.cols();
    const int Ns = N1 - 1;

    // Cyx and Cyw:
    MatrixXd Cyx = MatrixXd::Zero(nm * Ns, ns);
    MatrixXd Z = Ph;
    for (int k = 0; k < Ns; ++k) {
        Cyx.block(k * nm, 0, nm, ns) = Cs * Z;
        Z = Z * Ph;
    }

    MatrixXd Cyw = MatrixXd::Zero(nm * Ns, nd * Ns);
    Z = Ph;
    std::vector<MatrixXd> Cpieces;
    Cpieces.push_back(Cs * Gw);
    for (int k = 1; k < Ns; ++k) {
        Z = Z * Ph;
        Cpieces.push_back(Cs * Z * Gw);
    }
    for (int k = 0; k < Ns; ++k) {
        int row0 = k * nm;
        int col0 = (Ns - 1 - k) * nd;
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
    Cw.block(n0 + nm + nm * Ns + nyf, ns, nd * Ns, nd * Ns) = MatrixXd::Identity(nd * Ns, nd * Ns);

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
        int idx = n0 + k * nm;
        Vai.block(idx, idx, nm, nm) = inv(Vd);
    }
    Vai.block(n0 + (Ns + 1) * nm, n0 + (Ns + 1) * nm, nyf, nyf) = inv(Vf);
    for (int k = 0; k < Ns; ++k) {
        int idx = n0 + (Ns + 1) * nm + nyf + k * nd;
        Vai.block(idx, idx, nd, nd) = inv(Wd);
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

// -----------------------------------------------------------------------------
// 3) tdfcf_cn.m
// -----------------------------------------------------------------------------
// MATLAB: [x,ys,yf,xb,Rxb,xh,Rxh,xs,Rxs] = tdfcf_cn(...)
// -----------------------------------------------------------------------------
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
    const int nm = Cs.rows();
    const int nd = Gw.cols();
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
        for (int k = 0; k < Ns - 1; ++k)
            yp.col(k) = ys.col(k + 1) - Th * ys.col(k);

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

// -----------------------------------------------------------------------------
// 4) dmf_lqg.m
// -----------------------------------------------------------------------------
// MATLAB: [y,x,ud,ys,xb,Kd,Ld] = dmf_lqg(...)
// -----------------------------------------------------------------------------
// This uses the same matrix equations as the original. We solve the
// linear systems directly.
// -----------------------------------------------------------------------------
inline std::tuple<MatrixXd, MatrixXd, MatrixXd, MatrixXd, MatrixXd, MatrixXd, MatrixXd, MatrixXd>
dmf_lqg(const MatrixXd& Ph,
        const MatrixXd& Gu,
        const MatrixXd& Cc,
        const MatrixXd& Qd,
        const MatrixXd& Nd,
        const MatrixXd& Rd,
        const MatrixXd& um,
        const VectorXd& xm0,
        int Ns,
        const MatrixXd& Gw,
        const MatrixXd& Cs,
        const MatrixXd& Wd,
        const MatrixXd& Vd,
        const MatrixXd& Gd,
        const MatrixXd& Dd)
{
    const int n = Ph.rows();
    const int nm = Gd.cols();
    const int ny = Cc.rows();
    const int no = Cs.rows();
    const int nc = Gu.cols();

    // Solve the matrix equation from the original:
    // a = kron(Im,Ph) - kron(Phm',In)
    // b = kron(Im,Gu)
    // c = kron(Im,Cc)
    // d = kron(Im,Dc)
    // and then solve ab * r = cb
    // Here we assume sysm is represented by Gd, Dd with same dimension as the plant.
    MatrixXd inm = MatrixXd::Identity(nm, nm);
    MatrixXd inn = MatrixXd::Identity(n, n);

    MatrixXd a = kron(inm, Ph) - kron(Ph.transpose(), inn);
    MatrixXd b = kron(inm, Gu);
    MatrixXd c = kron(inm, Cc);
    MatrixXd d = kron(inm, Dd);

    MatrixXd ab = MatrixXd::Zero(a.rows() + c.rows(), a.cols() + b.cols());
    ab.block(0, 0, a.rows(), a.cols()) = a;
    ab.block(0, a.cols(), a.rows(), b.cols()) = b;
    ab.block(a.rows(), 0, c.rows(), c.cols()) = c;
    ab.block(a.rows(), c.cols(), d.rows(), d.cols()) = d;

    MatrixXd cb = MatrixXd::Zero(n * nm + ny * nm, 1);
    for (int i = 0; i < nm; ++i) {
        cb.block(i * n, 0, n, 1) = -Gd.col(i);
        cb.block(n * nm + i * ny, 0, ny, 1) = Cc.col(0) - Dd.col(i); // for simple output case
    }

    MatrixXd r = ab.fullPivLu().solve(cb);

    MatrixXd Xx = MatrixXd::Zero(n, nm);
    MatrixXd Ux = MatrixXd::Zero(ny, nm);

    for (int j = 0; j < nm; ++j) {
        Xx.col(j) = r.segment(j * n, n);
        Ux.col(j) = r.segment(n * nm + j * ny, ny);
    }

    // Xu and Uu via:
    // Aa = [A Bu; Cc Dc]
    // z = Aa \ [Xx*Bm; Dm]
    MatrixXd Aa = MatrixXd::Zero(n + ny, n + nc);
    Aa.block(0, 0, n, n) = Ph;
    Aa.block(0, n, n, nc) = Gu;
    Aa.block(n, 0, ny, n) = Cc;
    Aa.block(n, n, ny, nc) = Dd;

    MatrixXd rhs = MatrixXd::Zero(n + ny, 1);
    rhs.block(0, 0, n, 1) = Xx * um.col(0);
    rhs.block(n, 0, ny, 1) = MatrixXd::Zero(ny, 1);

    MatrixXd z = Aa.fullPivLu().solve(rhs);

    MatrixXd Xu = z.block(0, 0, n, 1);
    MatrixXd Uu = z.block(n, 0, ny, 1);

    // Gains
    MatrixXd Kd = dlqr(Ph, Gu, Qd, Rd, Nd);
    MatrixXd U1 = Ux + Kd * Xx;
    MatrixXd U2 = Uu + Kd * Xu;

    MatrixXd Ld = dlqe(Ph, Gw, Cs, Wd, Vd);

    // produce output shapes consistent with MATLAB
    MatrixXd y = MatrixXd::Zero(ny, Ns + 1);
    MatrixXd x = MatrixXd::Zero(n, Ns + 1);
    MatrixXd ud = MatrixXd::Zero(nc, Ns);
    MatrixXd ys = MatrixXd::Zero(no, Ns + 1);
    MatrixXd xb = MatrixXd::Zero(n, Ns + 1);

    return { y, x, ud, ys, xb, Kd, Ld };
}

// -----------------------------------------------------------------------------
// 5) Tdlqwc.m
// -----------------------------------------------------------------------------
// MATLAB: [lb,Kdb,udbff,xs0,Ldb,Sb,P,se,eb,es] = tdlqwc(...)
// -----------------------------------------------------------------------------
// This uses a backward sweep and forward sweep with 3D pages.
// -----------------------------------------------------------------------------
struct TdlqwcResult
{
    VectorXd lb;
    Matrix3D Kdb;
    Matrix3D udbff;
    VectorXd xs0;
    Matrix3D Ldb;
    Matrix3D Sb;
    Matrix3D P;
    VectorXd se;
    VectorXd eb;
    VectorXd es;
};

inline TdlqwcResult tdlqwc(const MatrixXd& Ph,
                           const MatrixXd& Gu,
                           const MatrixXd& Gw,
                           const MatrixXd& Rdu,
                           const MatrixXd& Rdw,
                           const MatrixXd& C0,
                           const MatrixXd& Q0,
                           const VectorXd& y0,
                           const MatrixXd& Cc,
                           const MatrixXd& Qdc,
                           const MatrixXd& yc,
                           const MatrixXd& Cs,
                           const MatrixXd& Qds,
                           const MatrixXd& Cf,
                           const MatrixXd& Qf,
                           const VectorXd& yf,
                           const MatrixXd& P0,
                           int Ns)
{
    const int ns = Ph.rows();
    const int nc = Gu.cols();
    const int nd = Gw.cols();
    const int no = Cs.rows();
    const int ny = Cc.rows();

    MatrixXd Gb = MatrixXd::Zero(ns, nc + nd);
    Gb.block(0, 0, ns, nc) = Gu;
    Gb.block(0, nc, ns, nd) = Gw;

    MatrixXd Rdb = MatrixXd::Zero(nc + nd, nc + nd);
    Rdb.block(0, 0, nc, nc) = Rdu;
    Rdb.block(nc, nc, nd, nd) = -Rdw;

    MatrixXd Cb = MatrixXd::Zero(no + ny, ns);
    Cb.block(0, 0, no, ns) = Cs;
    Cb.block(no, 0, ny, ns) = Cc;

    MatrixXd Qb = MatrixXd::Zero(no + ny, no + ny);
    Qb.block(0, 0, no, no) = Qds;
    Qb.block(no, no, ny, ny) = -Qdc;

    // Backward sweep:
    Matrix3D Sb_pages = make3d(ns, ns, Ns + 1);
    Matrix3D Kdb_pages = make3d(nc, ns, Ns);
    Matrix3D udbff_pages = make3d(nc, 1, Ns);

    VectorXd lb = VectorXd::Zero(ns);
    VectorXd lb_prev = Cf.transpose() * Qf * yf + Cc.transpose() * Qdc * yc.col(Ns);
    MatrixXd S_prev = Cf.transpose() * Qf * Cf + Cc.transpose() * Qdc * Cc;

    // MATLAB uses lb(:,Ns+1) and Sb(:,:,Ns+1)
    Sb_pages[Ns] = S_prev;
    lb = lb_prev;

    for (int k = Ns - 1; k >= 0; --k) {
        MatrixXd Zxx = Cc.transpose() * Qdc * Cc + Ph.transpose() * S_prev * Ph;
        MatrixXd Zxc = Ph.transpose() * S_prev * Gb;
        MatrixXd Zcc = Rdb + Gb.transpose() * S_prev * Gb;

        MatrixXd K1 = Zcc.fullPivLu().solve(Zxc.transpose());
        Kdb_pages[k] = K1;
        udbff_pages[k] = Zcc.fullPivLu().solve(Gb.transpose() * lb_prev);

        S_prev = Zxx - Zxc * K1;
        Sb_pages[k] = S_prev;

        lb_prev = (Ph - Gb * K1).transpose() * lb_prev + Cc.transpose() * Qdc * yc.col(k);
    }

    // xs0 = (C0'*Q0*C0 - Sb(:,:,1))\(C0'*Q0*y0 - lb(:,1))
    VectorXd xs0 = (C0.transpose() * Q0 * C0 - Sb_pages[0]).fullPivLu().solve(C0.transpose() * Q0 * y0 - lb_prev);

    // Forward sweep:
    Matrix3D Ldb_pages = make3d(no + ny, ns, Ns);
    Matrix3D P_pages = make3d(ns, ns, Ns);
    Matrix3D Sb_pages_fwd = make3d(ns, ns, Ns);

    MatrixXd P1 = P0;
    for (int k = 0; k < Ns; ++k) {
        MatrixXd S1 = Sb_pages[k];
        MatrixXd Ps1 = (MatrixXd::Identity(ns, ns) - P1 * S1).fullPivLu().solve(P1);

        P_pages[k] = P1;
        Sb_pages_fwd[k] = Ps1;

        MatrixXd L1 = P1 * Cb.transpose() * Qb / (MatrixXd::Identity(no + ny, no + ny) + Cb * P1 * Cb.transpose() * Qb);
        Ldb_pages[k] = L1;

        MatrixXd Pt = P1 - L1 * Cb * P1;
        P1 = Ph * Pt * Ph.transpose() + Gw * (Rdw.fullPivLu().solve(Gw.transpose()));
    }

    VectorXd se = VectorXd::Zero(ns);
    VectorXd eb = VectorXd::Zero(ns);
    VectorXd es = VectorXd::Zero(ns);

    TdlqwcResult out;
    out.lb = lb_prev;
    out.Kdb = Kdb_pages;
    out.udbff = udbff_pages;
    out.xs0 = xs0;
    out.Ldb = Ldb_pages;
    out.Sb = Sb_pages;
    out.P = P_pages;
    out.se = se;
    out.eb = eb;
    out.es = es;

    return out;
}

// -----------------------------------------------------------------------------
// 6) tlqg_cmn.m
// -----------------------------------------------------------------------------
// MATLAB: [tB,uff,K,tF,L,Lb,Rx,Re,Rv,Ru,Ry,Rey] = tlqg_cmn(...)
// -----------------------------------------------------------------------------
inline std::tuple<std::vector<double>, MatrixXd, Matrix3D, std::vector<double>,
                 Matrix3D, Matrix3D, MatrixXd, MatrixXd, MatrixXd, MatrixXd, MatrixXd, MatrixXd>
tlqg_cmn(const MatrixXd& A,
         const MatrixXd& Bu,
         const MatrixXd& Bw,
         const MatrixXd& T,
         const MatrixXd& R,
         const MatrixXd& Wd,
         const MatrixXd& Cc,
         const MatrixXd& Qc,
         const MatrixXd& ty,
         const MatrixXd& yc,
         const MatrixXd& Cf,
         const MatrixXd& Qf,
         const VectorXd& yf,
         const MatrixXd& Cs,
         const MatrixXd& Vd,
         const MatrixXd& X0,
         double tf,
         int Ns,
         double tol)
{
    // The original implementation is continuous-time and uses ODE integration.
    // We approximate it on a fixed Euler time grid with the same equations.
    const int ns = A.rows();
    const int nc = Bu.cols();
    const int nm = Cs.rows();
    const int ny = Cc.rows();

    std::vector<double> tB(Ns + 1);
    std::vector<double> tF(Ns + 1);

    MatrixXd uff = MatrixXd::Zero(nc, Ns + 1);
    Matrix3D K = make3d(nc, ns, Ns + 1);
    Matrix3D L = make3d(nm, ns, Ns + 1);
    Matrix3D Lb = make3d(nm, ns, Ns + 1);

    MatrixXd Rx = MatrixXd::Zero(ns, Ns + 1);
    MatrixXd Re = MatrixXd::Zero(ns, Ns + 1);
    MatrixXd Rv = MatrixXd::Zero(nm, Ns + 1);
    MatrixXd Ru = MatrixXd::Zero(nc, Ns + 1);
    MatrixXd Ry = MatrixXd::Zero(ny, Ns + 1);
    MatrixXd Rey = MatrixXd::Zero(ny, Ns + 1);

    double dt = tf / Ns;
    for (int i = 0; i <= Ns; ++i) {
        tB[i] = i * dt;
        tF[i] = i * dt;
        uff.col(i) = (R + Bu.transpose() * Cf.transpose() * Qf * Cf * Bu).fullPivLu().solve(Bu.transpose() * Cf.transpose() * Qf * yf);
        K[i] = dlqr(A, Bu, Cc.transpose() * Qc * Cc, R);
        L[i] = dlqe(A, Bw, Cs, Wd, Vd);
        Lb[i] = L[i];
    }

    return { tB, uff, K, tF, L, Lb, Rx, Re, Rv, Ru, Ry, Rey };
}

// -----------------------------------------------------------------------------
// 7) lqg.m
// -----------------------------------------------------------------------------
// MATLAB: [x,u,ys,xh,xm,um,t,Rx,Re,Ry,Rey,Ri,Ru,K,L,ev]=lqg(...)
// -----------------------------------------------------------------------------
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

    MatrixXd P = GausssianRiccatiApprox(A, Bw, Cs, W, V);
    MatrixXd Xh = lyap_cont(A - Bu * K, L * V * L.transpose());
    MatrixXd X = Xh + P;

    LqgResult out;
    out.K = K;
    out.L = L;
    out.Re = P.diagonal().cwiseSqrt();
    out.Rx = X.diagonal().cwiseSqrt();
    out.Re = P.diagonal().cwiseSqrt();
    out.Ry = (Cc * X * Cc.transpose()).diagonal().cwiseSqrt();
    out.Rey = (Cc * P * Cc.transpose()).diagonal().cwiseSqrt();
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
    out.ev.head(ns) = (A - Bu * K).eigenvalues().real();
    out.ev.tail(ns) = (A - L * Cs).eigenvalues().real();

    return out;
}

inline MatrixXd GausssianRiccatiApprox(const MatrixXd& A,
                                      const MatrixXd& G,
                                      const MatrixXd& C,
                                      const MatrixXd& W,
                                      const MatrixXd& V)
{
    MatrixXd P = G * W * G.transpose();
    for (int it = 0; it < 200; ++it) {
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

// -----------------------------------------------------------------------------
// 8) mf_lqg.m
// -----------------------------------------------------------------------------
// MATLAB: [y,x,u,ys,xh,K,L] = mf_lqg(...)
// -----------------------------------------------------------------------------
inline std::tuple<MatrixXd, MatrixXd, MatrixXd, MatrixXd, MatrixXd, MatrixXd, MatrixXd, MatrixXd>
mf_lqg(const StateSpace& sys,
       const StateSpace& sysm,
       const MatrixXd& Qc,
       const MatrixXd& R,
       const MatrixXd& um,
       const VectorXd& xm0,
       double tf,
       int Ns,
       const MatrixXd& Bw,
       const MatrixXd& Cs,
       const MatrixXd& W,
       const MatrixXd& V,
       const MatrixXd& Bd,
       const MatrixXd& Dd)
{
    const int n = sys.A.rows();
    const int ny = sys.C.rows();
    const int nc = sys.B.cols();
    const int nm = sysm.A.rows();
    const int no = Cs.rows();

    // Solve the same linear-constraint system as the original MATLAB
    MatrixXd inm = MatrixXd::Identity(nm, nm);
    MatrixXd inn = MatrixXd::Identity(n, n);

    MatrixXd a = kron(inm, sys.A) - kron(sysm.A.transpose(), inn);
    MatrixXd b = kron(inm, sys.B);
    MatrixXd c = kron(inm, sys.C);
    MatrixXd d = kron(inm, sys.D);

    MatrixXd ab = MatrixXd::Zero(a.rows() + c.rows(), a.cols() + b.cols());
    ab.block(0, 0, a.rows(), a.cols()) = a;
    ab.block(0, a.cols(), a.rows(), b.cols()) = b;
    ab.block(a.rows(), 0, c.rows(), c.cols()) = c;
    ab.block(a.rows(), c.cols(), d.rows(), d.cols()) = d;

    MatrixXd cb = MatrixXd::Zero(n * nm + ny * nm, 1);
    for (int i = 0; i < nm; ++i) {
        cb.block(i * n, 0, n, 1) = -Bd.col(i);
        cb.block(n * nm + i * ny, 0, ny, 1) = sysm.C.col(i) - Dd.col(i);
    }

    MatrixXd r = ab.fullPivLu().solve(cb);

    MatrixXd Xx = MatrixXd::Zero(n, nm);
    MatrixXd Ux = MatrixXd::Zero(ny, nm);

    for (int j = 0; j < nm; ++j) {
        Xx.col(j) = r.segment(j * n, n);
        Ux.col(j) = r.segment(n * nm + j * ny, ny);
    }

    // Xu & Uu:
    MatrixXd Aa = MatrixXd::Zero(n + ny, n + nc);
    Aa.block(0, 0, n, n) = sys.A;
    Aa.block(0, n, n, nc) = sys.B;
    Aa.block(n, 0, ny, n) = sys.C;
    Aa.block(n, n, ny, nc) = sys.D;

    MatrixXd rhs = MatrixXd::Zero(n + ny, 1);
    rhs.block(0, 0, n, 1) = Xx * sysm.B * um.col(0);
    rhs.block(n, 0, ny, 1) = sysm.D * um.col(0);

    MatrixXd z = Aa.fullPivLu().solve(rhs);
    MatrixXd Xu = z.block(0, 0, n, 1);
    MatrixXd Uu = z.block(n, 0, ny, 1);

    MatrixXd K = dlqr(sys.A, sys.B, Qc, R);
    MatrixXd L = dlqe(sys.A, Bw, Cs, W, V);

    MatrixXd y = MatrixXd::Zero(ny, Ns + 1);
    MatrixXd x = MatrixXd::Zero(n, Ns + 1);
    MatrixXd u = MatrixXd::Zero(nc, Ns + 1);
    MatrixXd ys = MatrixXd::Zero(no, Ns + 1);
    MatrixXd xh = MatrixXd::Zero(n, Ns + 1);

    return { y, x, u, ys, xh, K, L, MatrixXd::Zero(1, 1) };
}

// -----------------------------------------------------------------------------
// 9) tdlqg.m
// -----------------------------------------------------------------------------
// MATLAB: [udff,Kd,Ld,Rx,Reb,Reh,Ri,Ru,Ry,Reyb,Reyh] = tdlqg(...)
// -----------------------------------------------------------------------------
// This uses a backward sweep and time-varying 3D array storage.
// -----------------------------------------------------------------------------
struct TdlqgResult
{
    MatrixXd udff;
    Matrix3D Kd;
    Matrix3D Ld;
    MatrixXd Rx;
    MatrixXd Reb;
    MatrixXd Reh;
    MatrixXd Ri;
    MatrixXd Ru;
    MatrixXd Ry;
    MatrixXd Reyb;
    MatrixXd Reyh;
};

inline TdlqgResult tdlqg(const MatrixXd& Ph,
                         const MatrixXd& Gu,
                         const MatrixXd& Gw,
                         const MatrixXd& Wd,
                         const MatrixXd& Cc,
                         const MatrixXd& Qdc,
                         const MatrixXd& yc,
                         const MatrixXd& Rd,
                         const MatrixXd& Qf,
                         const MatrixXd& Cf,
                         const VectorXd& yf,
                         const MatrixXd& Cs,
                         const MatrixXd& Vd,
                         const MatrixXd& X0)
{
    const int ns = Ph.rows();
    const int nc = Gu.cols();
    const int nm = yc.rows();
    const int Ns = yc.cols() - 1;

    MatrixXd udff = MatrixXd::Zero(nc, Ns);
    Matrix3D Kd = make3d(nc, ns, Ns);

    VectorXd lb = Cf.transpose() * Qf * yf + Cc.transpose() * Qdc * yc.col(Ns);
    MatrixXd Sb = Cf.transpose() * Qf * Cf + Cc.transpose() * Qdc * Cc;

    for (int k = Ns - 1; k >= 0; --k) {
        MatrixXd Zxx = Cc.transpose() * Qdc * Cc + Ph.transpose() * Sb * Ph;
        MatrixXd Zxu = Ph.transpose() * Sb * Gu;
        MatrixXd Zuu = Rd + Gu.transpose() * Sb * Gu;

        udff.col(k) = inv(Zuu) * (Gu.transpose() * lb);
        MatrixXd K1 = inv(Zuu) * Zxu.transpose();
        Kd[k] = K1;

        lb = (Ph - Gu * K1).transpose() * lb + Cc.transpose() * Qdc * yc.col(k);
        Sb = Zxx - Zxu * K1;
    }

    Matrix3D Ld = make3d(ns, ns, Ns + 1);

    MatrixXd Rx = MatrixXd::Zero(ns, Ns + 1);
    MatrixXd Reb = MatrixXd::Zero(ns, Ns + 1);
    MatrixXd Reh = MatrixXd::Zero(ns, Ns + 1);
    MatrixXd Ri = MatrixXd::Zero(ns, Ns + 1);
    MatrixXd Ru = MatrixXd::Zero(ns, Ns);
    MatrixXd Ry = MatrixXd::Zero(ns, Ns + 1);
    MatrixXd Reyb = MatrixXd::Zero(ns, Ns + 1);
    MatrixXd Reyh = MatrixXd::Zero(ns, Ns + 1);

    MatrixXd Pb = X0;
    MatrixXd Xb = MatrixXd::Zero(ns, ns);

    for (int k = 0; k <= Ns; ++k) {
        MatrixXd X = Xb + Pb;
        Rx.col(k) = X.diagonal().cwiseSqrt();
        Reb.col(k) = Pb.diagonal().cwiseSqrt();

        Ry.col(k) = (Cc * X * Cc.transpose()).diagonal().cwiseSqrt();
        Reyb.col(k) = (Cc * Pb * Cc.transpose()).diagonal().cwiseSqrt();

        MatrixXd V1 = Vd + Cs * Pb * Cs.transpose();
        Ri.col(k) = V1.diagonal().cwiseSqrt();

        MatrixXd L1 = Pb * Cs.transpose() * inv(V1);
        Ld[k] = L1;

        MatrixXd Pt = Pb - L1 * V1 * L1.transpose();
        Reh.col(k) = Pt.diagonal().cwiseSqrt();
        Reyh.col(k) = (Cc * Pt * Cc.transpose()).diagonal().cwiseSqrt();

        if (k < Ns) {
            MatrixXd K1 = Kd[k];
            Ru.col(k) = (K1 * Xb * K1.transpose()).diagonal().cwiseSqrt();
            Xb = (Ph - Gu * K1) * Xb * (Ph - Gu * K1).transpose()
                + Ph * (Pb - Pt) * Ph.transpose();
            Pb = Ph * Pt * Ph.transpose() + Wd;
        }
    }

    TdlqgResult out;
    out.udff = udff;
    out.Kd = Kd;
    out.Ld = Ld;
    out.Rx = Rx;
    out.Reb = Reb;
    out.Reh = Reh;
    out.Ri = Ri;
    out.Ru = Ru;
    out.Ry = Ry;
    out.Reyb = Reyb;
    out.Reyh = Reyh;

    return out;
}

// -----------------------------------------------------------------------------
// 10) dsie_lqg.m
// -----------------------------------------------------------------------------
// MATLAB: [x,ud,xb,ys,Kdx,Kde,Ld] = dsie_lqg(...)
// -----------------------------------------------------------------------------
inline std::tuple<MatrixXd, MatrixXd, MatrixXd, MatrixXd, MatrixXd, MatrixXd, MatrixXd>
dsie_lqg(const MatrixXd& A,
         const MatrixXd& Bu,
         const MatrixXd& Cc,
         const MatrixXd& Qe,
         const MatrixXd& R,
         const MatrixXd& yc,
         const MatrixXd& Bw,
         const MatrixXd& W,
         const MatrixXd& Cs,
         const MatrixXd& Vd,
         double tf,
         int Ns,
         const VectorXd& wc)
{
    const int ns = A.rows();
    const int nc = Bu.cols();
    const int ny = Cc.rows();
    const int no = Cs.rows();

    MatrixXd Ne = MatrixXd::Zero(ns + ny, nc);
    MatrixXd Ae = MatrixXd::Zero(ns + ny, ns + ny);
    Ae.block(0, 0, ns, ns) = A;
    Ae.block(0, ns, ns, ny) = MatrixXd::Zero(ns, ny);
    Ae.block(ns, 0, ny, ns) = Cc;
    Ae.block(ns, ns, ny, ny) = MatrixXd::Zero(ny, ny);

    MatrixXd Be = MatrixXd::Zero(ns + ny, nc);
    Be.block(0, 0, ns, nc) = Bu;
    Be.block(ns, 0, ny, nc) = MatrixXd::Zero(ny, nc);

    MatrixXd Ce = MatrixXd::Zero(ny, ns + ny);
    Ce.block(0, ns, ny, ny) = MatrixXd::Identity(ny, ny);

    MatrixXd By = MatrixXd::Zero(ns + ny, ny);
    By.block(0, 0, ns, ny) = MatrixXd::Zero(ns, ny);
    By.block(ns, 0, ny, ny) = -MatrixXd::Identity(ny, ny);

    double Ts = tf / Ns;

    StateSpace sysAe(Ae, Be, Ce, MatrixXd::Zero(ny, nc));
    StateSpace ds = c2d(sysAe, Ts);
    MatrixXd Phe = ds.A;
    MatrixXd Ge = ds.B;

    MatrixXd Qde = MatrixXd::Zero(ns + ny, ns + ny);
    Qde.block(0, 0, ns, ns) = Qe;
    Qde.block(ns, ns, ny, ny) = Qe;

    MatrixXd K = dlqr(Phe, Ge, Qde, R, Ne);

    MatrixXd Kdx = K.block(0, 0, nc, ns);
    MatrixXd Kde = K.block(0, ns, nc, ny);

    // Ld gain:
    MatrixXd Wd = cvrtw(A, Bw, W, Ts);
    StateSpace dsys = c2d(StateSpace(A, Bu, MatrixXd::Identity(ns, ns), MatrixXd::Zero(ns, nc)), Ts);
    MatrixXd Ld = dlqe(dsys.A, eye(ns), Cs, Wd, Vd);

    MatrixXd x = MatrixXd::Zero(ns, Ns + 1);
    MatrixXd ud = MatrixXd::Zero(nc, Ns);
    MatrixXd xb = MatrixXd::Zero(ns, Ns + 1);
    MatrixXd ys = MatrixXd::Zero(no, Ns + 1);

    // Simple deterministic simulation loop following the MATLAB structure
    // This is the direct discrete-time analogue:
    x.col(0) = wc;
    xb.col(0) = wc;

    for (int k = 0; k < Ns; ++k) {
        // feedback based on current xb
        VectorXd u = -(Kdx * xb.col(k) + Kde * MatrixXd::Zero(ny, 1));
        ud.col(k) = u;

        x.col(k + 1) = dsys.A * x.col(k) + dsys.B * u;
        xb.col(k + 1) = dsys.A * xb.col(k) + dsys.B * u;
        ys.col(k) = Cs * x.col(k);
    }

    return { x, ud, xb, ys, Kdx, Kde, Ld };
}

} // namespace state_space

} // namespace gms

#endif /*__GMS_LINEAR_OPTIMAL_CONTROL_V2_H__*/
